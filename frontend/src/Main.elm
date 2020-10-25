port module Main exposing (..)

import Array exposing (Array)
import Browser exposing (Document)
import Browser.Events exposing (onKeyDown)
import Debug
import Elastic exposing (..)
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (onClick, onInput)
import Http exposing (..)
import Info exposing (introduction)
import Json.Decode as Decode
import Keyboard.Event exposing (KeyboardEvent, considerKeyboardEvent, decodeKeyboardEvent)
import Nav exposing (navbar)
import Process
import Regex
import Task exposing (..)


port consoleLog : String -> Cmd msg


type alias SearchSuggestions =
    Array String


type alias ActiveSuggestion =
    Int


type alias Model =
    { searchString : String
    , searchSuggestions : SearchSuggestions
    , activeSuggestion : Maybe Int
    }


init : ( Model, Cmd Msg )
init =
    ( { searchString = ""
      , searchSuggestions = Array.empty
      , activeSuggestion = Nothing
      }
    , Cmd.none
    )



---- UPDATE ----


type Key
    = ArrowUp
    | ArrowDown


type Msg
    = SearchUpdate String
    | Search
    | SearchResponse (Result Http.Error String)
    | GetSearchSuggestions String
    | GotSearchSuggestions (Result Http.Error SearchSuggestions)
    | KeyPressed Key
    | SuggestionSelected Int


maybeSerialize result =
    -- Just used for console logging during development
    case result of
        Err value ->
            Debug.toString value

        Ok value ->
            serialize value


delay : Float -> msg -> Cmd msg
delay time msg =
    Process.sleep time
        |> Task.andThen (always <| Task.succeed msg)
        |> Task.perform identity


decodeSearchSuggestions : Decode.Decoder SearchSuggestions
decodeSearchSuggestions =
    Decode.field "suggestions" (Decode.array Decode.string)


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    case msg of
        SearchUpdate value ->
            if value == "" then
                -- Don't wait to clear the suggestions if there is nothing
                -- in the search box. Having search suggestions with an empty
                -- searchString causes issues for the highlightMatchingText function
                ( { model
                    | searchString = value
                    , activeSuggestion = Nothing
                    , searchSuggestions = Array.empty
                  }
                , Cmd.none
                )

            else
                ( { model
                    | searchString = value
                    , activeSuggestion = Nothing
                  }
                , delay 1000 (GetSearchSuggestions value)
                )

        Search ->
            let
                consoleMsg =
                    parse model.searchString |> maybeSerialize

                serverQuery =
                    get { url = "/search/" ++ consoleMsg, expect = Http.expectString SearchResponse }
            in
            ( model, Cmd.batch [ consoleLog consoleMsg, serverQuery ] )

        SearchResponse result ->
            ( model, Debug.toString result |> consoleLog )

        GetSearchSuggestions value ->
            -- This is some debounce on the search string to prevent spamming ElasticSearch with queries
            if model.searchString == value then
                ( model
                , get
                    { url = "/search_as_you_type/" ++ value
                    , expect = Http.expectJson GotSearchSuggestions decodeSearchSuggestions
                    }
                )

            else
                ( model, Cmd.none )

        GotSearchSuggestions result ->
            ( { model | searchSuggestions = Result.withDefault Array.empty result }, Debug.toString result |> consoleLog )

        KeyPressed key ->
            let
                activeSuggestion =
                    if model.searchSuggestions /= Array.empty then
                        Just
                            (modBy (Array.length model.searchSuggestions)
                                (case model.activeSuggestion of
                                    Nothing ->
                                        case key of
                                            ArrowUp ->
                                                Array.length model.searchSuggestions

                                            ArrowDown ->
                                                0

                                    Just value ->
                                        case key of
                                            ArrowUp ->
                                                value - 1

                                            ArrowDown ->
                                                value + 1
                                )
                            )

                    else
                        Nothing
            in
            ( { model | activeSuggestion = activeSuggestion }, consoleLog (Debug.toString activeSuggestion) )

        SuggestionSelected value ->
            let
                searchString =
                    case Array.get value model.searchSuggestions of
                        Just string ->
                            string

                        Nothing ->
                            model.searchString
            in
            ( { model | searchSuggestions = Array.empty, searchString = searchString }, consoleLog searchString )



---- VIEW ----


highlightMatchingText : String -> String -> List (Html Msg)
highlightMatchingText searchString suggestion =
    -- Note splitting a string removes the split string eg.> split "a" "James" = ["J", "mes"]
    if searchString /= "" then
        -- Determine location of matches (Case insensitive!)
        String.indexes (String.toLower searchString) (String.toLower suggestion)
            -- Get List of matching characters
            |> List.map (\idx -> String.slice idx (idx + (String.length searchString)) suggestion)
            -- Consecutively split the suggestion on the matching strings
            |> List.foldl (\str -> List.concatMap (\innerStr -> (String.split str innerStr))) [suggestion]
            -- Convert each split string to NON-bold text
            |> List.map (\t -> p [ class "d-inline" ] [ text t ])
            -- Insert BOLD text of matching strings
            |> List.intersperse (p [ class "d-inline font-weight-bold" ] [ text searchString ])

    else
        []


viewSuggestions : String -> SearchSuggestions -> Maybe ActiveSuggestion -> Html Msg
viewSuggestions searchString suggestions active =
    let
        selector =
            case active of
                Nothing ->
                    \_ -> ""

                Just value ->
                    \idx ->
                        if idx == value then
                            "active"

                        else
                            ""
    in
    ul [ class "list-group" ]
        (List.indexedMap
            (\idx suggestion ->
                li
                    [ String.join " "
                        [ "list-group-item list-group-item-action"
                        , selector idx
                        ]
                        |> class
                    , onClick (SuggestionSelected idx)
                    ]
                    (highlightMatchingText searchString suggestion)
            )
            (Array.toList suggestions)
        )


view : Model -> Document Msg
view model =
    { title = "Digital Expression Explorer 2"
    , body =
        [ navbar
        , div [ class "container my-5 mx-5 mx-auto" ]
            [ input
                [ onInput SearchUpdate
                , attribute "aria-label" "Search"
                , class "form-control form-control-lg"
                , placeholder "Human epilepisy | SRP070529"
                , type_ "search"
                , value model.searchString
                ]
                []
            , viewSuggestions model.searchString model.searchSuggestions model.activeSuggestion
            -- Alternate ^^^'viewSuggestions' func with no highlighting useful to debug
            --, ul [] (Array.toList (Array.map (\str -> li [][text str]) model.searchSuggestions))
            , div [ class "d-flex justify-content-center" ]
                [ button
                    [ onClick Search
                    , class "btn btn-lg btn-outline-success my-5"
                    , type_ "button"
                    ]
                    [ text "Search" ]
                ]
            ]
        , introduction
        ]
    }


toKey : KeyboardEvent -> Maybe Msg
toKey keyboardEvent =
    -- IF Nothing is returned the update
    -- function is never called which simplifies
    -- downstream logic
    case keyboardEvent.key of
        Just value ->
            if value == "ArrowUp" then
                Just (KeyPressed ArrowUp)

            else if value == "ArrowDown" then
                Just (KeyPressed ArrowDown)

            else
                Nothing

        _ ->
            Nothing


subscriptions : Model -> Sub Msg
subscriptions model =
    onKeyDown (considerKeyboardEvent toKey)



---- PROGRAM ----


main : Program () Model Msg
main =
    Browser.document
        { view = view
        , init = \_ -> init
        , update = update
        , subscriptions = subscriptions
        }

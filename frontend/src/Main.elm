port module Main exposing (..)

import Array
import Browser exposing (Document)
import Browser.Events exposing (onKeyDown)
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (onClick)
import Info exposing (introduction)
import KeyBoardHelpers exposing (enterKey)
import Keyboard.Event exposing (KeyboardEvent, considerKeyboardEvent)
import Nav exposing (navbar)
import SearchBar
import SearchBarViews exposing (..)


port consoleLog : String -> Cmd msg


setPage : Model -> Page -> Model
setPage model page =
    { model | page = page }


type Page
    = Home
    | SearchResults


type alias Model =
    { searchBar : SearchBar.Model
    , page : Page
    }


init : ( Model, Cmd Msg )
init =
    ( { searchBar = SearchBar.init
      , page = Home
      }
    , Cmd.none
    )



---- UPDATE ----


type Msg
    = GotSearchBarMsg SearchBar.Msg
    | Search -- Message of this type will be sent to the SearchBar module
    | EnterKey


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        fromSearchBar =
            \( mdl, cmd ) -> ( { model | searchBar = mdl }, Cmd.map GotSearchBarMsg cmd )
    in
    case msg of
        GotSearchBarMsg message ->
            SearchBar.update message model.searchBar
                |> fromSearchBar

        Search ->
            let
                ( mdl, cmd ) =
                    SearchBar.update SearchBar.searchMsg model.searchBar
                        |> fromSearchBar
            in
            ( setPage mdl SearchResults, cmd )

        EnterKey ->
            case ( Array.isEmpty model.searchBar.searchSuggestions, model.searchBar.activeSuggestion ) of
                ( False, Just value ) ->
                    -- Selecting a suggestion with enter takes
                    -- precedence over searching with enter
                    SearchBar.update (SearchBar.suggestionSelectedMsg value) model.searchBar
                        |> fromSearchBar

                ( _, _ ) ->
                    -- search with enter
                    SearchBar.update SearchBar.searchMsg model.searchBar
                        |> fromSearchBar
                        |> (\( mdl, cmd ) -> ( { mdl | page = SearchResults }, cmd ))



---- VIEW ----


searchButton : Html Msg
searchButton =
    div [ class "d-flex justify-content-center" ]
        [ button
            [ onClick Search
            , class "btn btn-lg btn-outline-success my-5"
            , type_ "button"
            ]
            [ text "Search" ]
        ]


pageLayout : List (Html Msg) -> List (Html Msg)
pageLayout content =
    [ navbar
    , div [ class "container my-5 mx-5 mx-auto" ] content
    , introduction
    ]


pageView : Model -> List (Html Msg)
pageView model =
    let
        fromSearchBar =
            Html.map GotSearchBarMsg
    in
    pageLayout <|
        case model.page of
            Home ->
                [ fromSearchBar (viewLargeSearchBar model.searchBar), searchButton ]

            SearchResults ->
                [ fromSearchBar (viewSearchResults model.searchBar.searchResults) ]


view : Model -> Document Msg
view model =
    { title = "Digital Expression Explorer 2"
    , body = pageView model
    }


subscriptions : Model -> Sub Msg
subscriptions model =
    Sub.batch
        [ Sub.map GotSearchBarMsg (SearchBar.subscriptions model.searchBar)
        , onKeyDown (considerKeyboardEvent (enterKey EnterKey))
        ]



---- PROGRAM ----


main : Program () Model Msg
main =
    Browser.document
        { view = view
        , init = \_ -> init
        , update = update
        , subscriptions = subscriptions
        }
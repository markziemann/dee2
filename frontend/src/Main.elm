port module Main exposing (..)

import Array
import Browser exposing (Document)
import Browser.Events exposing (onKeyDown)
import Browser.Navigation as Nav
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (onClick)
import Info exposing (introduction)
import KeyBoardHelpers exposing (arrowDown, arrowUp, enterKey, try)
import Keyboard.Event exposing (KeyboardEvent, considerKeyboardEvent)
import MainTypes exposing (..)
import MainViews
import Nav exposing (navbar)
import SearchBar
import SearchBarTypes
import SearchBarViews exposing (..)
import Table
import Url
import Url.Parser as UrlP exposing ((</>), (<?>))
import Url.Parser.Query as Query


port consoleLog : String -> Cmd msg


homeSlug =
    "/"


searchResultsSlug =
    "SearchResults"


routeParser : UrlP.Parser (Route -> a) a
routeParser =
    UrlP.oneOf
        [ UrlP.map HomeRoute (UrlP.s homeSlug)
        , UrlP.map SearchResultsRoute (UrlP.s searchResultsSlug <?> Query.string "q")
        ]


homePage =
    HomePage HomeRoute


searchResultsPage : Maybe String -> Page
searchResultsPage maybeString =
    SearchResultsPage (SearchResultsRoute maybeString)


determinePage : Url.Url -> Page
determinePage url =
    case UrlP.parse routeParser url of
        Just page ->
            case page of
                HomeRoute ->
                    homePage

                SearchResultsRoute maybeString ->
                    searchResultsPage maybeString

        Nothing ->
            homePage


init : () -> Url.Url -> Nav.Key -> ( Model, Cmd Msg )
init flags url navKey =
    ( { navKey = navKey
      , url = url
      , searchBar = SearchBar.init
      , page = determinePage url
      , searchHits = Nothing
      , searchResultRows = Nothing
      , resultsTableState = Table.initialSort "id"
      , resultsTableQuery = ""
      }
    , Cmd.none
    )



---- UPDATE ----


requestSearch model fromSearchBar =
    SearchBar.update SearchBar.searchMsg model.searchBar
        |> fromSearchBar
        |> (\( mdl, cmd ) ->
                ( mdl
                , Cmd.batch
                    [ [ searchResultsSlug, model.searchBar.searchString ]
                        |> String.join "/?q="
                        |> Nav.pushUrl model.navKey
                    , cmd
                    ]
                )
           )


updateSearchData : Model -> SearchBarTypes.SearchResults -> Model
updateSearchData model searchResults =
    { model | searchHits = Just searchResults.hits, searchResultRows = Just searchResults.rows }


maybeUpdate : (Model -> a -> Model) -> Maybe a -> Model -> Model
maybeUpdate updateFunc maybe model =
    case maybe of
        Just value ->
            updateFunc model value

        Nothing ->
            model


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        fromSearchBar =
            \( mdl, cmd, searchResults ) ->
                ( { model | searchBar = mdl } |> maybeUpdate updateSearchData searchResults
                , Cmd.map GotSearchBarMsg cmd
                )
    in
    case msg of
        GotSearchBarMsg message ->
            SearchBar.update message model.searchBar
                |> fromSearchBar

        Search ->
            requestSearch model fromSearchBar

        EnterKey ->
            case ( Array.isEmpty model.searchBar.searchSuggestions, model.searchBar.activeSuggestion ) of
                ( False, Just value ) ->
                    -- Selecting a suggestion with enter takes
                    -- precedence over searching with enter
                    SearchBar.update (SearchBar.suggestionSelectedMsg value) model.searchBar
                        |> fromSearchBar

                ( _, _ ) ->
                    requestSearch model fromSearchBar

        ResultClicked result ->
            let
                searchResultRows =
                    Maybe.map (Array.set result.id { result | selected = not result.selected })
                        model.searchResultRows
            in
            ( { model | searchResultRows = searchResultRows }, Cmd.none )

        SetResultsTableQuery resultsTableQuery ->
            ( { model | resultsTableQuery = resultsTableQuery }
            , Cmd.none
            )

        SetResultsTableState resultsTableState ->
            ( { model | resultsTableState = resultsTableState }
            , Cmd.none
            )

        LinkClicked urlRequest ->
            case urlRequest of
                Browser.Internal url ->
                    ( model, Nav.pushUrl model.navKey (Url.toString url) )

                Browser.External href ->
                    ( model, Nav.load href )

        UrlChanged url ->
            ( { model | url = url, page = determinePage url }
            , Cmd.none
            )



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
            HomePage pageData ->
                [ fromSearchBar (viewLargeSearchBar model.searchBar), searchButton ]

            SearchResultsPage pageData ->
                [ MainViews.viewSearchResults model ]


view : Model -> Document Msg
view model =
    { title = "Digital Expression Explorer 2"
    , body = pageView model
    }


subscriptions : Model -> Sub Msg
subscriptions model =
    case model.page of
        HomePage route ->
            Sub.batch
                [ Sub.map GotSearchBarMsg SearchBar.subscriptions
                , onKeyDown (considerKeyboardEvent (enterKey EnterKey))
                ]

        SearchResultsPage route ->
            Sub.none



---- PROGRAM ----


main : Program () Model Msg
main =
    Browser.application
        { view = view
        , init = init
        , update = update
        , subscriptions = subscriptions
        , onUrlRequest = LinkClicked
        , onUrlChange = UrlChanged
        }

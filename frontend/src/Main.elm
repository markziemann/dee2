port module Main exposing (..)

import Browser exposing (Document)
import Html exposing (..)
import Html.Events exposing (onClick)
import Html.Attributes exposing (..)
import Info exposing (introduction)
import Nav exposing (navbar)
import SearchBar
import SearchBarViews exposing (..)


port consoleLog : String -> Cmd msg


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


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        fromSearchBar = (\( a, b ) -> ( { model | searchBar = a }, Cmd.map GotSearchBarMsg b ))
    in
    case msg of
        GotSearchBarMsg message ->
            SearchBar.update message model.searchBar
                |> fromSearchBar
        Search ->
            SearchBar.update SearchBar.searchMsg model.searchBar
                |> fromSearchBar



---- VIEW ----

searchButton: Html Msg
searchButton =
    div [ class "d-flex justify-content-center" ]
                [ button
                    [ onClick Search
                    , class "btn btn-lg btn-outline-success my-5"
                    , type_ "button"
                    ]
                    [ text "Search" ]
                ]

pageLayout: List (Html Msg) -> List (Html Msg)
pageLayout content =
    [ navbar
    , div [ class "container my-5 mx-5 mx-auto" ] content
    , introduction
    ]



pageView : Model -> List (Html Msg)
pageView model =
    let
        fromSearchBar = Html.map GotSearchBarMsg
    in
    pageLayout <| case model.page of
        Home ->
            [fromSearchBar (viewLargeSearchBar model.searchBar), searchButton]
        SearchResults ->
            [fromSearchBar (viewSearchResults model.searchBar.searchResults)]


view : Model -> Document Msg
view model =
    { title = "Digital Expression Explorer 2"
    , body = pageView model
    }


subscriptions : Model -> Sub Msg
subscriptions model =
    Sub.map GotSearchBarMsg (SearchBar.subscriptions model.searchBar)



---- PROGRAM ----


main : Program () Model Msg
main =
    Browser.document
        { view = view
        , init = \_ -> init
        , update = update
        , subscriptions = subscriptions
        }
